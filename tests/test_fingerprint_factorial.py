"""Controls for label semantics, preprocessing and paired factorial diagnostics."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import torch
from tests.test_fingerprint_review import data_fixture
from meta_acdc.structure.receptor_identity_audit import seqres_sequences,same_construct_conflicts,is_subsequence
from meta_acdc.models.fingerprint_factorial import assemble,epoch_order,ARMS
from scripts import run_fingerprint_factorial as runner
import scripts.build_fingerprint_dataset as builder

class IdentityEvidenceChecks(unittest.TestCase):
    def test_seqres_rejects_incomplete_or_unknown(self):
        def line(chain,n,names):
            chars=[' ']*80;chars[:6]='SEQRES';chars[11]=chain;chars[13:17]=f'{n:4d}';chars[19:19+len(names)]=names
            return ''.join(chars)
        text='\n'.join([line('D',3,'ALA ARG GLY'),line('E',4,'ALA ARG GLY'),line('I',2,'ALA UNK')])
        self.assertEqual(seqres_sequences(text),{'D':'ARG'})
        self.assertTrue(is_subsequence('AG','ARG'))
        self.assertFalse(is_subsequence('AT','ARG'))
    def test_equivalence_requires_both_complete_chains(self):
        a=dict(eligible_construct_evidence=True,observed_sequence_id='first',seqres_alpha='A',seqres_beta='B')
        b=dict(a,observed_sequence_id='second')
        self.assertEqual(len(same_construct_conflicts([a,b])),1)
        self.assertEqual(same_construct_conflicts([a,dict(b,seqres_beta='C')]),[])
        self.assertEqual(same_construct_conflicts([a,dict(b,eligible_construct_evidence=False)]),[])
    def test_builder_requires_explicit_diagnostic_mode(self):
        with patch('sys.argv',['build','--inventory','unused.pt']),contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit) as error:
            builder.main()
        self.assertEqual(error.exception.code,2)

class FactorialChecks(unittest.TestCase):
    def test_identity_only_restores_named_block(self):
        d=data_fixture();a=assemble(d['graphs'],d['records'],[0,1]);b=assemble(d['graphs'],d['records'],[0,1],True)
        raw=torch.cat([torch.as_tensor(d['graphs'][i][0]) for i in [0,1]])
        torch.testing.assert_close(b[0][:,4:24],raw[:,4:24])
        torch.testing.assert_close(a[0][:,:4],b[0][:,:4]);torch.testing.assert_close(a[0][:,24:],b[0][:,24:])
        for i in [1,2,3,4,5]:torch.testing.assert_close(a[i],b[i])
    def test_batch_order_preserves_population_and_steps(self):
        d=data_fixture();idx=list(range(8))
        for grouped in [False,True]:
            x=epoch_order(d['records'],idx,random.Random(8),grouped)
            self.assertEqual(sorted(x),idx);self.assertEqual(len(set(x)),len(idx))
            self.assertFalse(set(x)&set(range(8,12)))
            self.assertEqual(x,epoch_order(d['records'],idx,random.Random(8),grouped))
    def test_numeric_geometry_mask_keeps_topology(self):
        d=data_fixture();a=assemble(d['graphs'],d['records'],[0]);b=assemble(d['graphs'],d['records'],[0],intervention='numeric_geometry_zero')
        self.assertTrue(torch.equal(a[2],b[2]));self.assertTrue(torch.equal(a[0][:,:25],b[0][:,:25]))
        self.assertEqual(float(b[1].abs().sum()+b[3].abs().sum()+b[0][:,25].abs().sum()),0)
    def test_tiny_factorial_replay_and_no_biological_promotion(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);d=data_fixture();d['source_indices']=list(range(12));torch.save(d,root/'dev.pt')
            p=dict(protocol_id='synthetic-unit-test',scope='synthetic observed-key diagnostic only',dataset=str(root/'dev.pt'),
                   frozen_sha256={str(root/'dev.pt'):hashlib.sha256((root/'dev.pt').read_bytes()).hexdigest()},
                   fit_indices=list(range(8)),validation_indices=list(range(8,12)),quarantined_source_indices=[99],
                   known_same_pairs=[[0,2],[2,0]],seeds=[7,11,17],interventions=['identity_zero','numeric_geometry_zero','peptide_pool_zero'],
                   recipe=dict(threads=1,emb_dim=8,lr=3e-4,lr_min=3e-6,weight_decay=1e-4,epochs=2,warmup_epochs=1,batch_size=4,temperature=0.07))
            path=root/'protocol.json';path.write_text(json.dumps(p));path.with_suffix('.sha256').write_text(hashlib.sha256(path.read_bytes()).hexdigest())
            with patch('sys.argv',['run','--protocol',str(path),'--out-dir',str(root/'run')]),contextlib.redirect_stdout(io.StringIO()):runner.main()
            result=json.loads((root/'run/evaluation.json').read_text());status=json.loads((root/'run/status.json').read_text())
            self.assertEqual(status['state'],'complete');self.assertEqual(status['n_main_models'],12)
            self.assertFalse(result['biological_receptor_candidate_promotable']);self.assertTrue(all(result['same_seed_initial_weights_equal'].values()))
            self.assertEqual(set(result['aggregate']),set(ARMS));self.assertIn('interaction',result['factorial_effects']['validation']['raw_hard'])

if __name__=='__main__':unittest.main()
