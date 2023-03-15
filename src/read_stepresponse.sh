#!/bin/sh
# Extract data from Google Drive using your favorite method. MATLAB's webread may be an alternative method. Google Drive may ask for spam confirmation for files bigger than 40 MB. Sometimes throws ERROR: cannot verify docs.google.com's certificate, but still works.

# Or, download manually at:
# https://drive.google.com/file/d/1-Uxd52iJMb0VudF8Q5t5vrd7Vnqm6UO1/view?usp=sharing

# Download ttdsummary_global_13july2011.mat from Google Drive
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1-Uxd52iJMb0VudF8Q5t5vrd7Vnqm6UO1' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1-Uxd52iJMb0VudF8Q5t5vrd7Vnqm6UO1" -O ttdsummary_global_13july2011.mat && rm -rf /tmp/cookies.txt
