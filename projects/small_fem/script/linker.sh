#!/bin/bash

echo $1 | sed -e 's|/.*/lib|-l|' -e 's|.so||'
