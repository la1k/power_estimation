"""
Continuously poll the specified rotctld instance for angles,
and print timestamp and angle if angles change.

Usage:

python3 rotctld_angle_printer.py ROTCTLD_HOST[:ROTCTLD_PORT] > angle_file

Example: `python3 rotctld_angle_printer.py localhost` will connect to the
rotctld instance on localhost:4533 (default rotctld port) and print the rotor
angles to stdout.
"""

import socket
import time
import datetime
import sys

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit()

#parse command line options
host_port_args = sys.argv[1].split(':') #assume argument 1 on form host:port
rotctld_host = host_port_args[0]
rotctld_port = 4533
if len(host_port_args) > 1:
    rotctld_port = int(host_port_args[1])

#polling time
poll_time_ms = 10

#connect to rotctld
rotctl = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
rotctl.connect((rotctld_host, rotctld_port))

prev_azimuth = float('NaN')
prev_elevation = float('NaN')

print("timestamp\tazimuth\televation")

while True:
    time.sleep(poll_time_ms*1.0/(1000.0))

    #get current angle
    rotctl.send(b'p\n')
    azimuth, elevation = rotctl.recv(1024).decode('ascii').splitlines()

    #print if it differs from previous angle
    if (azimuth != prev_azimuth) or (elevation != prev_elevation):
        print(datetime.datetime.now().isoformat() + "\t" + str(azimuth) + '\t' + str(elevation), flush=True)

    prev_azimuth = azimuth
    prev_elevation = elevation


