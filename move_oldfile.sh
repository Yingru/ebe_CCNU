#!/bin/bash
# rm files older than 60 minutes in results directory
find results/ -type f -mmin +60 -exec rm -f '{}' +
