#!/usr/bin/env python3
import urllib.request
import xml.etree.ElementTree as ET

sras = [
    'SRR28369364','SRR28369365','SRR28369367','SRR28369368',
    'SRR28369369','SRR28369370','SRR28369371','SRR28369373',
    'SRR28369375','SRR28369376','SRR28369377','SRR28369378',
    'SRR28369380','SRR28369381','SRR28369382','SRR28369385','SRR28369386'
]

url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={','.join(sras)}&rettype=runinfo&retmode=xml"
req = urllib.request.urlopen(url)
tree = ET.fromstring(req.read())

print("SampleID\tGroup\tCohort\tBodySite")
for row in tree.findall('.//Row'):
    run = row.find('Run').text
    sname = row.find('SampleName').text
    if sname.startswith('UC'):
        group = 'Ulcerative_Colitis'
        cohort = 'UC (IBD)'
    elif sname.startswith('N'):
        group = 'Healthy_Control'
        cohort = 'Healthy'
    else:
        group = 'Cohort'
        cohort = sname
    print(f"{run}\t{group}\t{cohort}\tGut_Stool")
