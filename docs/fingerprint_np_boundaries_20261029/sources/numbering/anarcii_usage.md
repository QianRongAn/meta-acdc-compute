<br>

When the model is instantiated it can be done minimally:
```python
from anarcii import Anarcii

model = Anarcii()

seqs = ["EIVMTQSPDTLSVSPGERATLSCRASESISSNLAWYQQKPGQVPRLLIYGASTRATGVPARFTGSGSGTEFTLTISSLQSEDFAVYYCQQYNNRLPYTFGQGTKLEIKRTVAAP",
        "DIVMTQSRDTLSVSPGERATLSCRSSESISSNLEWYQQKPGQVPRLLIYGISTRATGVPQRFTGSGSGTQFTLTISSLQSQDFQVYYCQQYNNRLPYTFGQGTKLEIKRTVQQP"]

results = model.number(seqs)
```

<br>

Or the arguments can be used to modify the model being run, inference speed and output data formats.
```python
from anarcii import Anarcii

model = Anarcii(seq_type="tcr", 
                mode="speed", 
                batch_size=128, 
                cpu=True, 
                ncpu=32, 
                verbose=True, 
                max_seqs_len=10_000_000)

seqs = "lots_of_seqs.fasta"

results = model.number(seqs)
```

<br>

# Model Arguments

- **`seq_type`** *(str, default="antibody")*:   
  Defines the type of sequence being numbered. One of `"antibody"`, `"tcr"`, `"unknown"` or `"shark"`. This in combination with mode (below) determines the model being run.

- **`mode`** *(str, default="accuracy")*:  
  In combination with **`seq_type`* this determines the model. One of `"accuracy"` or `"speed"`.

- **`batch_size`** *(int, default=8)*:  
  Number of sequences processed at a time. See section [recommended batch sizes](https://github.com/ALGW71/ANARCII-DEV/wiki/FAQs#recommended-batch-sizes).

- **`cpu`** *(bool, default=False)*:  
  If set to `True`, explicitly forces the use of CPU computation instead of GPU, regardless of availability.

- **`ncpu`** *(int, default=-1)*:  
  Number of CPU cores to use. This does impact speed when running on a GPU (queuing of the dataloader), therefore it is recommended to increase this as much as possible to maximise speed when running on large numbers of sequences.

- **`verbose`** *(bool, default=False)*:  
  Provides more info on the progress during model processing. Useful for monitoring, especially when running on many sequences.

- **`max_seqs_len`** *(int, default=102400)*:  
  Defined the threshold for entering batch mode. If the input fasta file or list of sequences is greater than 102,400 then the model will enter batch mode to save RAM (writing outputs to a **`.msgpack`** file which can be accessed in a memory efficient way). Outputs must be accessed by running **`model.to_csv("myseqs.csv")`** etc or interacting with the **`.msgpack`** file. This can be increased if you are not limited by RAM and want the model to return an object in python. See [MessagePack](https://github.com/oxpig/ANARCII/wiki/Working-with-MessagePack).

- **`return_logits`** *(bool, default=False)*:  
  Returns the logit values associated with each position to allow users to inspect relative differences in model confidence. Please note that these are heavily context dependent, therefore it is challenging to compare the same values across different sequence types/content.

- **`low_score_warn`** *(bool, default=False)*:
Warns users when a sequence score falls outside the defined window from known sequence sets. See [Using ANARCII as a screening tool](https://github.com/oxpig/ANARCII/wiki/Using-ANARCII-as-a-screening-tool). Please note this feature is not available in batch mode.



