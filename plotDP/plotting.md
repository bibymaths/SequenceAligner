
## Visualizing with `gnuplot`

### Why use it?

* DP matrix and traceback files are large, often matching the product of the sequence lengths, making manual inspection impractical.
* Visual representation helps interpret alignment decisions and traceback paths more easily.
* `gnuplot` is a robust and flexible tool for generating high-quality heatmaps and matrix plots from plain text data.
 
- Install `gnuplot` if not already installed. For example, on Ubuntu: 

```bash
sudo apt-get install gnuplot
```  
 
or on Fedora:  

```bash 
sudo dnf install gnuplot
``` 

or on macOS: 

```bash
brew install gnuplot
``` 
 
- Run the script to generate plot from traceback/DP files  
 
```bash 
chmod +x plotDP.sh 
./plotDP.sh <lcs_traceback_file> <global_dp_file> <local_dp_file> 
```

