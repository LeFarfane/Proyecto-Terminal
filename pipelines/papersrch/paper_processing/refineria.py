import pandas as pd
import os
# Set up the path and read the file into the script.
path = r"D:\AAA\2. Estudios\Ingeniería en Biotecnología\Proyecto Terminal\Pipelines\data\papers\papers.jsonl"
papers = pd.read_json(path, lines=True)

#Creating safe copy
work = papers.copy(deep=True)

#let's find out if we have duplicated papers.
# check for same PMID

'''
duplicated_pmid = work[work["PMID"].duplicated(keep=False)]
print(duplicated_pmid.shape)
print(duplicated_pmid.head())
'''
#remove all the duplicated ones, by default it keeps only the first instance.
work = work.drop_duplicates(subset=["PMID"])
print("Numer of duplicated PMIDs:")
print(work["PMID"].duplicated().sum())
#print(work.shape)
#print(work.loc[0:5, "Title"])

#first we will normalize all titles so it's easier to filter through them.
work["Title"] = work["Title"].str.lower()
work["Abstract"] = work["Abstract"].str.lower()

#Now lets remove the undesired paper.
#We will do that by using a list of keywords that of no interest to us.

unwanted_keywords = [
    'dental',
    'cancer',
    'mouse',
    'rat',
    'mice',
    'drugs',
    'drug',
    'assessment',
    'ct',
    'surgery',
    'magnetic',
    'imagining',
    'mri',
    'vitamin',
    'herbal',
    'covid',
    'nutritional',
    'tumor',
    'acupuncture',
    'pharmacologic',
    'pharmacology',
    'postoperative',
    'coffee',
    'coping',
    'management',
    'diet',
    'oral',
    'heart',
    'stroke',
    'carcinogenesis',
    'oil',
    'psychological',
    'safety',
    'rhinitis',
    'consumption',
    'dosing',
    'coronavirus',
    'radiation',
    'vaccination',
    'violence',
    'telemedicine',
    ]


#now we detect which articles have the unwanted words. and will add a tag into another column to identify them.
def contains_bad(title, keywords):
    return any(keywords in title for keywords in keywords)

work["unwanted"] = work["Title"].apply(lambda x: contains_bad(x, unwanted_keywords))

#then we split it
dumpster = work[work["unwanted"]].copy()
work = work[~work["unwanted"]].copy()

print("Kept:", len(work), " → Dumpster:", len(dumpster))

#custumize the path
dumpster_path = r"D:\AAA\2. Estudios\Ingeniería en Biotecnología\Proyecto Terminal\Pipelines\data\papers\dumpster\dumpster.csv"

def save_to_dumpster(df, path):
    if not os.path.exists(path):
        # create new file with header if the file doesn't exists
        df.to_csv(path, index=False, mode="w", header=True)
        print(f"Created new dumpster file with {len(df)} rows → {path}")
        return
    else:
        # append rows without writing header to an existing file THIS SHOULD BE APPEND BUT CHANGED CUZ IT WAS ACCUMULATING AND WAS TOO LAZY TO FIX IT
        df.to_csv(path, index=False, mode="w", header=False)
        print(f"Appended {len(df)} rows to existing dumpster → {path}")
        return

#The file with unwanted papers is created
save_to_dumpster(dumpster, dumpster_path)



print(f"This is the info for work: {work.shape}")
print(f"This is the info for dumpster: {dumpster.shape}")
#print(f"This is the info for dumpster: {dumpster.shape}")

"""
    At this point we already separated most of the unrelated papers, now we will continue to do the opposite
what will do is to create a list of the most useful keywords and create a dataset with the most useful papers
"""
useful_keywords = [
    #,
    'microarn',
    'prediction',
    'ibd',
    'markers',
    'biomarkers',
    'artificial',
    'intelligence',
    ]


#now we detect which articles have the useful words. and will add a tag into another column to identify them.
def contains_useful(row, keywords):
    text = f"{row["Title"]} {row["Abstract"]}"
    return any(k in text for k in keywords)

work["wanted"] = work.apply(lambda row: contains_useful(row, useful_keywords), axis=1)

#then we split it
wanted = work[work["wanted"]].copy()
work = work[~work["wanted"]].copy()

print("Kept:", len(work), " → Useful:", len(wanted))

#custumize the path
useful_path = r"D:\AAA\2. Estudios\Ingeniería en Biotecnología\Proyecto Terminal\Pipelines\data\papers\dumpster\useful.csv"
wanted["Title"] = wanted["Title"].str.title()
save_to_dumpster(wanted, useful_path)

work_path = r"D:\AAA\2. Estudios\Ingeniería en Biotecnología\Proyecto Terminal\Pipelines\data\papers\dumpster\work.csv"
work["Title"] = work["Title"].str.title()
save_to_dumpster(work, work_path)