from spyprot.fetchRcsbChainInfo import getIdenticalChains, getSimilarChains, getUniqChains

a = getIdenticalChains("2jlo", chain="A").get()
print(a)
a = getSimilarChains("2jlo",chain="A",identity=90).get()
print(a)
a = getUniqChains("2jlo").get()
print(a)

