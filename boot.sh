#!/bin/bash

source venv/bin/activate
ping
gunicorn -b :5000 --timeout 300 --access-logfile - --error-logfile - app:app
