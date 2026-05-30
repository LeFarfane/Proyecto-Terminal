import pandas as pd
from pathlib import Path

path = Path("mirna_lit_seq.csv")

base =  {
    "mirna" : "string",
    "fold_value" : "float64",
    "source" : "string",
    "disease" : "string",
    "effect" : "string",
}

df = pd.DataFrame({c: pd.Series(dtype=t) for c, t in base.items() })
print(df.dtypes)

path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(path, index=False)

