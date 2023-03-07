# Smith Waterman Algorithm

Implements the Smith-Waterman algorithm for local sequence alignment. 

## Installation
Use the command 
```bash
pip install git+https://github.com/intermezzo672/smith-waterman.git
```

## Usage
```bash
from smith_waterman_intermezzo672 import hw1

# returns output txt file with sequences, score matrix, and alignment results 
hw1.runSW(input_file, score_file, openGap=-2, extGap=-1)
```
View an example output in ```/samples/sample-output1.txt```  
Input file should contain two sequences on two different lines   
View an example input in ```/samples/sample-input1.txt```  
View an example score file in ```/samples/blosum62.txt```

## Repository Command-Line Usage
```bash
python hw1.py -i <input_file> -s <score_file>
```

Example (assuming cd-ed into src directory):
```bash
python hw1.py -i ../samples/sample-input1.txt -s ../samples/blosum62.txt
```
