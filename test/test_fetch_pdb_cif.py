from fetchRcsbChainInfo import PdbFile, MMCIFfile

p = PdbFile("/tmp", "1k36", "A").download()
p = PdbFile("/tmp","1k36", "A", 'CA').download()
#p = MMCIFfile("/tmp", "1k36", "A").download()

# Download and extract chain 5
p = PdbFile("/tmp", "6p5i", "5").download()
# Download and extract only C3' atoms from chain 5
p = PdbFile("/tmp", "6p5i", "5", 'C3\'').download()

#m = MMCIFfile("/tmp", "1k36", "A").download()
# Download and extract chain 5
#m = MMCIFfile("/tmp", "6p5i", "5").download()
# Download and extract only C3' atoms from chain 5
m = MMCIFfile("/tmp", "6p5i", "5", 'C3\'').download()