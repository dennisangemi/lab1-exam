#!/bin/bash

# this script will download data from google sheets

curl -kL "https://docs.google.com/spreadsheets/d/e/2PACX-1vTP1WVVpoOLCSUYdQ6ebtHJMoCCNo5Gx25k4Vb6tE2QDtj0HR_G1gAHD4G_lWC0wBVKA5vBO-Xkoxuy/pub?gid=88278913&single=true&output=csv" > data/exp-data-1.csv
echo "exp-data-1.csv downloaded!"

curl -kL "https://docs.google.com/spreadsheets/d/e/2PACX-1vTP1WVVpoOLCSUYdQ6ebtHJMoCCNo5Gx25k4Vb6tE2QDtj0HR_G1gAHD4G_lWC0wBVKA5vBO-Xkoxuy/pub?gid=370513976&single=true&output=csv" > data/exp-data-2.csv
echo "exp-data-2.csv"