#!/bin/bash

sudo find . -type d -exec chmod 771 "{}" \;
sudo find . -type d -exec chmod g+s "{}" \;
sudo find . -type f -exec chmod 664 "{}" \;

