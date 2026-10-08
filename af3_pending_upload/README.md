# af3_pending_upload —— 待提交 AlphaFold3 的结构任务

来源仓库：`https://github.com/QianRongAn/meta-acdc-compute`（分支 `main`，commit `1b7f129`）
生成时间：2026-10-07

## 这个文件夹是什么

`af3_batch1/` + `af3_batch2/` 里**还没提交过 AF3** 的任务，已按批次拆好。每个 `.fasta` 是一个
**5 链复合物输入**（`alpha` / `beta` / `peptide` / `mhc` / `b2m`），可以直接整段复制粘贴进
AlphaFold Server 网页版（alphafoldserver.com）的输入框——一次一个 job。

## 数字（两套独立台账交叉验证，完全一致）

| 项 | 数量 | 口径 |
|---|---|---|
| batch1 + batch2 任务总数 | 392 | 50 + 342（不含 manifest） |
| 已提交（有 AF3 结果） | 36 | `data/raw/af3_instances/**` 与 `data/raw/af3_predictions/**` 里存在结果文件 |
| └ batch1 已提交 | 29 | 50 的 58% |
| └ batch2 已提交 | 7 | 342 的 2% |
| **已被 relabel 覆盖（不要重传）** | **9** | 见 `obsolete_relabeled_do_not_reupload.tsv` |
| **本次待传** | **347** | batch1 **12** + batch2 **335** |

验证方式：`data/processed/af3_qc.tsv`（质检台账，覆盖 batch1 29 个 / batch2 7 个）与根据
`data/raw/af3_instances` 结果文件统计出的集合**逐条一致**，说明「已传」的判定没有漏算或错算。

## 那 9 个为什么不算待传

`data/processed/af3_relabel_map.tsv` 记录了一次 ID 更正：部分任务的真实 PDB 不是原名对应的结构，
旧 ID 被改名成新 ID。这 9 个的**新 ID 已经提交并有结果**，所以旧 ID 的文件属于作废件：

| 旧 ID（batch 文件夹里的文件名） | 已提交为 |
|---|---|
| `1ao7_LGYGFVNYI` | `1qrn_LGYGFVNYI` / `1qse_LGYGFVNYI` |
| `1ao7_LLFGFPVYV` | `1qrn_LLFGFPVYV` / `1qse_LLFGFPVYV` |
| `1ao7_LLFGKPVYV` | `1qrn_LLFGKPVYV` |
| `1ao7_LLFGPVYV` | `1qrn_LLFGPVYV` / `1qse_LLFGPVYV` |
| `1g6r_EQYKFYSV` | `1tcr_EQYKFYSV` |
| `1g6r_GGAPWNPAMMI` | `1tcr_GGAPWNPAMMI` |
| `1g6r_QLSPFPFDL` | `1tcr_QLSPFPFDL` |
| `1mwa_SIYRYYGL` | `1tcr_SIYRYYGL` |
| `1oga_PKYVKQNTLKLAT` | `2vlj_PKYVKQNTLKLAT`（内容逐字节相同） |

> 注意 `1ao7_LGYGFVNYI.fasta` 与 `1qrn_LGYGFVNYI.fasta` 内容**不同**——前者是改名前的旧受体版本，
> 提交它等于交一次错结构。仓库自带的 `data/processed/af3_upload_todo.txt` 也是把它们排除在外的。

## 文件夹内容

```
af3_pending_upload/
├── af3_batch1/        12 个 .fasta（batch1 剩余全部）
├── af3_batch2/       335 个 .fasta
├── pending_upload_manifest.tsv  待传清单：batch / job_id / pdb / peptide / cdr3 / mhc_hint
├── already_uploaded.tsv         已提交的 36 个（核对用）
├── obsolete_relabeled_do_not_reupload.tsv  9 个作废件（不要传）
└── README.md
```

待传按 PDB 分组共 111 个结构，其中 `1qsf`、`2gj6`、`3d39`、`3d3v`、`3h9s`、`3kxf`、`3pwp`、
`3qfj`、`3uts`、`3utt`、`4ftv`、`4jrx`、`5c07`、`5c08`、`5c09`、`5c0a`、`5c0b`、`5c0c`、`5hyj`
各 8 个（共 152 个），是最大的一批。

## 提交建议

1. **按 PDB 分批**：同一 PDB 的多个肽共用同一套 TCR/MHC 链，一起提交便于对齐检查。
2. **AF3 网页版有配额**，347 个不可能一次做完；`pending_upload_manifest.tsv` 可直接按 `pdb` 排序取下一批。
3. **提交后一定抽查一个结构**：本仓库踩过坑——AF3 网页版曾把 5 条链拼成 1 条（707 残基），
   导致首批 17 个任务作废。检查方法：MHC 重链应约 270+ 残基、TCR 链 CA 数 ~440，不是单一长链。
4. 顺带留意 `mhc_hint` 为空的 job（如 `1mwa_*`、`2ckb_*` 等）——没有 MHC 限制提示，`peptide` 与
   晶体里那条肽长度不一致时更值得核对。

## 结果回收（入库）

AF3 跑完拿到 `folds_<时间戳>` 目录后，按仓库既有流程入库：

```bash
python src/meta_acdc/structure/import_af3_instances.py \
    --src <folds_目录> \
    --list data/processed/kn5_submission_list.tsv \
    --existing data/raw/af3_predictions \
    --out data/raw/af3_instances \
    --manifest data/processed/af3_instance_manifest.tsv

python src/meta_acdc/structure/af3_qc.py \
    --root data/raw/af3_instances \
    --out data/processed/af3_qc.tsv
```

`import_af3_instances.py` 会把每次独立提交放进单独的 instance 目录（B1、B2……），
目的是量化重提方差；跑完 `af3_qc.py` 后 `af3_qc.tsv` 会自动膨胀，`af3_upload_todo.txt` 也可以按
本 README 的口径重新生成。
