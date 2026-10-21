"""Run the pinned tracked test modules plus publication evidence tests."""
from pathlib import Path
import json,subprocess,sys,unittest
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'src'))
raw=subprocess.check_output(['git','-C',str(ROOT),'ls-tree','-r','--name-only','15eb5c60b6907b9ba4f472dad3da612d7366c559','tests'],text=True)
paths=[p for p in raw.splitlines() if p.startswith('tests/test_') and p.endswith('.py')]+['tests/test_preprint_evidence.py']
suite=unittest.TestSuite()
for p in paths:suite.addTests(unittest.defaultTestLoader.discover(str(ROOT/'tests'),pattern=Path(p).name))
result=unittest.TextTestRunner(verbosity=2).run(suite)
receipt=dict(base_test_commit='15eb5c60b6907b9ba4f472dad3da612d7366c559',modules=paths,tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),concurrent_uncommitted_test_modules_excluded=True,passed=result.wasSuccessful())
target=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'snapshot_test_receipt.json'
with target.open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
sys.exit(not result.wasSuccessful())
