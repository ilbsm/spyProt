import urllib.request

get_page = urllib.request.urlopen('https://www.uniprot.org/uniprot/Q01292.txt')
get_ver = get_page.readlines(100)

print("Currant Versions Are: ", get_ver[0])
