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

def get_host_port(hostport_str):
    """
    Convert string in format `host:port` to host and port.
    """
    host_port_args = hostport_str.split(':') #assume argument 1 on form host:port
    rotctld_host = host_port_args[0]
    rotctld_port = 4533
    if len(host_port_args) > 1:
        rotctld_port = int(host_port_args[1])
    return rotctld_host, rotctld_port

def rotctld_connect(rotctld_host, rotctld_port):
    """
    Connect to rotctld server.
    """
    rotctl = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    rotctl.connect((rotctld_host, rotctld_port))
    return rotctl

def rotctld_angle(rotctl):
    """
    Get azimuth,elevation from rotctld server connection.
    """
    rotctl.send(b'p\n')
    azimuth, elevation = rotctl.recv(1024).decode('ascii').splitlines()
    return azimuth, elevation

def print_angles_forever(rotctl, fileobj=sys.stdout):
    """
    Poll angles indefinitely from rotctl, output to file.
    """
    poll_time_ms = 10

    prev_azimuth = float('NaN')
    prev_elevation = float('NaN')

    #header
    fileobj.write("timestamp\tazimuth\televation\n")

    while True:
        time.sleep(poll_time_ms*1.0/(1000.0))

        #get current angle
        azimuth, elevation = rotctld_angle(rotctl)

        #print if it differs from previous angle
        if (azimuth != prev_azimuth) or (elevation != prev_elevation):
            fileobj.write(datetime.datetime.now().isoformat() + "\t" + str(azimuth) + '\t' + str(elevation) + "\n")
        fileobj.flush()

        prev_azimuth = azimuth
        prev_elevation = elevation

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit()

    #parse input arguments
    host, port = get_host_port(sys.argv[1])

    #connect to rotctld
    rotctl = rotctld_connect(host, port)

    #print angles to stdout
    print_angles_forever(rotctl)
