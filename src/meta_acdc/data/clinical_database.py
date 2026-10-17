"""ClinicalTrials.gov aggregate safety evidence, never peptide-binding labels."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .release import sha256


def normalize_trial(study,source):
    protocol=study['protocolSection'];ident=protocol['identificationModule'];nct=ident['nctId'];module=study.get('resultsSection',{}).get('adverseEventsModule',{});groups={x['id']:x for x in module.get('eventGroups',[])};rows=[]
    for category in ['seriousEvents','otherEvents']:
        for event in module.get(category,[]):
            for stat in event.get('stats',[]):
                if stat['groupId'] not in groups:raise ValueError('unknown trial event group')
                affected=stat.get('numAffected');at_risk=stat.get('numAtRisk')
                if affected is not None and (affected<0 or (at_risk is not None and affected>at_risk)):raise ValueError('invalid affected/at-risk count')
                rows.append({'nct':nct,'group_id':stat['groupId'],'group_title':groups[stat['groupId']]['title'],'category':category,'term':event['term'],'organ_system':event.get('organSystem'),'num_affected':affected,'num_at_risk':at_risk,'num_events':stat.get('numEvents'),'source_vocabulary':event.get('sourceVocabulary'),'assessment_type':event.get('assessmentType'),'time_frame':module.get('timeFrame'),'reporting_frequency_threshold':module.get('frequencyThreshold'),'source':source,'binding_label':None,'peptide':None,'causal_offtarget_label':None,'unit':'trial arm, not individual patient/receptor-peptide pair'})
    return {'nct':nct,'title':ident.get('officialTitle',ident['briefTitle']),'last_update':protocol.get('statusModule',{}).get('lastUpdatePostDateStruct'),'interventions':protocol.get('armsInterventionsModule',{}).get('interventions',[]),'event_groups':list(groups.values()),'description':module.get('description'),'time_frame':module.get('timeFrame'),'source':source,'individual_patient_data':False,'paired_receptor_sequence_available_in_event_rows':False,'offtarget_peptide_available_in_event_rows':False},rows


def molecular_label(row):
    raise ValueError('Aggregate clinical adverse events cannot label a receptor-peptide binding or causal off-target pair')


def build(snapshot,output):
    sources=json.loads((snapshot/'retrieval_manifest.json').read_text());tracked={x['file']:x for x in sources};trials=[];events=[]
    for p in sorted(snapshot.glob('clinical_NCT*.json')):
        item=tracked[p.name]
        if sha256(p)!=item['sha256']:raise ValueError('changed clinical source')
        trial,rows=normalize_trial(json.loads(p.read_text()),item);trials.append(trial);events.extend(rows)
    audit=[]
    for nct in ['NCT01350401','NCT01352286']:
        trial=next(x for x in trials if x['nct']==nct)
        audit.append({'nct':nct,'historical_reference_pmid':'23926201','historical_reference_construct':'MAGE-A3 a3a / HLA-A*01','historical_reference_url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6002776/','current_interventions':trial['interventions'],'status':'current_intervention_mismatch_with_historical_reference','action':'Do not join current NY-ESO-1 outcome rows to historical MAGE-A3/titin cases by NCT alone; version-specific construct adjudication required','original_primary_xml_capture':'unavailable: Europe PMC HTTP500; historical reference inspected through primary web text, not counted as a locally hashed primary XML'})
    output.mkdir(parents=True,exist_ok=False)
    report={'n_trials':len(trials),'n_arm_event_observations':len(events),'n_event_groups':sum(len(x['event_groups']) for x in trials),'clinical_model_validation_allowed':False,'uses':'retrospective aggregate safety and identity audit only','limitations':['No individual patient dose/exposure/follow-up to receptor-peptide mapping','All-cause mortality is not treatment-caused mortality','AE occurrence does not establish molecular off-target causality','Unreported events and absent rows do not mean zero risk','No synthetic safe negative for a peptide or TCR','Do not sum counts across overlapping trials or event categories','Anti-CD3 soluble TCR bispecific and TCR-T mechanisms must be stratified']}
    for name,value in [('trial_summaries.json',trials),('arm_event_observations.json',events),('construct_link_audit.json',audit),('summary.json',report)]: (output/name).write_text(json.dumps(value,indent=2)+'\n')
    (output/'manifest.json').write_text(json.dumps({'sources':[x for x in sources if x['file'].startswith('clinical_')],'files':{p.name:sha256(p) for p in output.iterdir()},'clinical_model_validation_allowed':False},indent=2)+'\n');return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(build(a.snapshot,a.out),indent=2))
if __name__=='__main__':main()
