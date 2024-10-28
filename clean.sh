#!/bin/bash

find . -type d -name "unused" -exec git rm -r --cached {} +