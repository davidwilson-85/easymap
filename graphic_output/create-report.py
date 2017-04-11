#!/usr/bin/python

# This program dinamically creates an HTML-formatted report that serves both for the
# command-line and client-server versions of Easymap


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-project_name', action="store", dest='project_name')
args = parser.parse_args()
project_name = args.project_name

# Define file that contains the report
report_file = open(project_name + '/3_workflow_output/report.html', 'w')

# ONLY NEEDED FOR COMMAND-LINE VERSION. Write HTML document header to file
#report_file.write('<!DOCTYPE html>\n<html>\n<title>Easymap</title>\n<meta charset="UTF-8">\n<body>\n')

# Write beginning of report content
report_file.write('<div id="report">\n')

report_file.write('<h3 style="font-family:arial;">Section</h3>\n')

# Use this path for command-line version
report_file.write('<img src="test-image.jpg" width="700"></img><br>\n')
# Use this path for client-server version
report_file.write('<img src="../' + project_name +'/3_workflow_output/test-image.jpg" width="700"></img><br>\n')
# Temporary. For testing only
report_file.write('<img src="../../../interface/test-image.jpg" width="700"></img><br>\n')

# Write end of report content
report_file.write('</div>\n')

# ONLY NEEDED FOR COMMAND-LINE VERSION. Write HTML closing tags
#report_file.write('</body>\n</html>')

report_file.close()
















