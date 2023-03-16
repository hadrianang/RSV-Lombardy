# RSV-Lombardy
Repository for model code and aggregated data for the manuscript "Reconstructing the impact of COVID-19 on the immunity gap and transmission of respiratory syncytial virus in Lombardy, Italy." 

```
RSV-Lombardy
├── ModelInputFiles
|   ├── InputData.rds               #Input data matrices (in CSV format below)
|   ├── Hospitalisations.csv        #Hospital discharges [1]
|   ├── Population.csv              #Lombardy Population (ISTAT [2])
|   ├── PosTest.csv                 #Number of positive RSV tests [3]
|   ├── RSV_att_ILI.csv             #Number of RSV-attributable ILI cases [3]
|   └── TotalTests.csv              #Number of RSV tests performed [3]
├── CatalyticModel.stan             #Stan model code
├── RunModelLocal.R                 #R code for fitting Stan model
└── README.md
```
### References
[1] Regione Lombardia. Weekly hospital discharges in Lombardy - Unpublished. 2022 

[2] Italian National Institutde of Statistics. Lombardia Population. Resid. Popul. 1st January. 2022. http://dati.istat.it/?lang=en (accessed April 28, 2022).

[3] Istituto Superiore di Sanità. Sistema di Sorveglianza Integrata dell’Influenza. InfluNet. https://w3.iss.it/site/RMI/influnet/pagine/stagioni.aspx (accessed Aug 9, 2022).
