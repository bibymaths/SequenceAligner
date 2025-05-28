  
# Split hg38 genome into individual chromosome files 

```bash 
zcat hg38.fa.gz | awk '
  /^>/ { 
    if (seq) {
      print seq > (out)
    }
    out = substr($1, 2) ".fa"
    print $0 > out
    seq = ""
  }
  /^[^>]/ {
    seq = seq $0
  }
  END {
    if (seq) print seq > out
  }
' 
``` 
 
# Move chromosome files into a directory
```bash 
mkdir -p chr
mv chr* chr/  
```