# Physics lab 1 exam
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Experimental data and report of physics lab 1 exam

2022-09-05

## Repository structure
```
# root
## 📂 data
### 📄 exp-data-1.csv
### 📄 exp-data-2.csv
## 📂 report
### 📄 report.tex
### 📄 report.pdf
## 📂 scripts
## 📄 datapackage.yaml
## 📄 README.md

root
├── 📂 data
│   ├── 📄 exp-data-1.csv
│   └── 📄 exp-data-2.csv
├── 📂 report
│   ├── 📄 report.tex
│   └── 📄 report.pdf
├── 📂 scripts
├── 📄 datapackage.yaml
└── 📄 README.md
```

## Data dictionary and metadata

Every file in `data` folder is described according to [frictionless data standards](https://frictionlessdata.io/standards/) in the [`datapackage.yaml`](https://github.com/indecis-it/data/blob/main/datapackage.yaml)

### [📄 file.csv](https://github.com/indecis-it/data/blob/main/data/sources.csv)

- Path: `data/`
- URL: https://raw.githubusercontent.com/indecis-it/data/main/data/sources.csv
- Delimiter: `,`
- Encoding: `UTF-8`

field | type | description | example
-- | -- | -- | --
id | integer | Source ID | 1
title | string | Source title | Programma elettorale PD

### Unit of measurement
Units of measurement (column `uom`) are represented by alphanumeric codes according to [UNECE Recommendation 20](https://datahub.io/core/unece-units-of-measure). For more information you can visit the web page linked above.

## How to access data
If you wanna use these data in your works, you can follow these easy instructions. You have to locate the raw URL of the file you are interested in and paste it inside the formulas. Remember to read the [license](#license) and cite the author

### MATLAB
```matlab
categories = readtable("https://raw.githubusercontent.com/indecis-it/data/main/data/categories.csv")
```

### R
```r
categories <- read.csv("https://raw.githubusercontent.com/indecis-it/data/main/data/categories.csv")
```

## License
Data, code and images are licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). <br> <br>
<a href="https://creativecommons.org/licenses/by/4.0/"><img src="https://mirrors.creativecommons.org/presskit/buttons/88x31/png/by.png" width="150"/></a>

## Author
**Dennis Angemi**

Physics student at University of Catania (DFA)

📫 How to reach me:
  - [Email](mailto:dennisangemi@gmail.com)
  - [Twitter](https://twitter.com/dennisangemi)
  - [Telegram](https://t.me/dennisangemi)
  - [Instagram](http://instagram.com/dennisangemi)
