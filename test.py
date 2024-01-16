import urllib.request

from spyprot import PdbMetaData, UniprotSearch

get_page = urllib.request.urlopen('https://www.uniprot.org/uniprot/Q01292.txt')
get_ver = get_page.readlines(100)

print("Currant Versions Are: ", get_ver[0])


print(PdbMetaData('2efv').get())

print(UniprotSearch(fields=['accession', 'xref_pdb'], query='xref:pdb-2efv OR xref:pdb-1j85').get())