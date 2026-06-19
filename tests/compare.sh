#!/usr/bin/env bash
diff -s <(sort "$1") <(sort "$2")
