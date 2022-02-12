#!/bin/bash

docker run --rm -v $(pwd):/usr/share/nginx/html:ro -p 8082:80 nginx
