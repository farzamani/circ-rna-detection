import os
import csv


def parse_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    result = {}
    for line in lines:
        if line.startswith("\tCommand being timed:"):
            result['command'] = line.split(": ")[1].strip()
        elif line.startswith("\tUser time (seconds):"):
            result['user_time'] = float(line.split(": ")[1])
        elif line.startswith("\tSystem time (seconds):"):
            result['system_time'] = float(line.split(": ")[1])
        elif line.startswith("\tPercent of CPU this job got:"):
            result['cpu_percent'] = int(line.split(": ")[1].replace('%', ''))
        elif line.startswith("\tElapsed (wall clock) time"):
            elapsed_time_str = line.split(": ")[1].strip()
            elapsed_time_parts = elapsed_time_str.split(':')
            if len(elapsed_time_parts) == 2:
                elapsed_time_parts = ['0'] + elapsed_time_parts
            elapsed_time_seconds = float(elapsed_time_parts[2])
            elapsed_time_seconds = int(elapsed_time_seconds)
            elapsed_time_str = f"{elapsed_time_parts[0]}:{elapsed_time_parts[1]}:{elapsed_time_seconds:02d}"
            result['elapsed_time'] = elapsed_time_str
        elif line.startswith("\tMaximum resident set size"):
            result['max_resident_set_size'] = int(line.split(": ")[1])
            
    return result


def append_to_csv(csv_file, data):
    header =    [  
                    'command', 
                    'user_time', 
                    'system_time', 
                    'elapsed_time', 
                    'cpu_percent', 
                    'max_resident_set_size'
                ]
    
    # Create file if it does not exist
    if not os.path.exists(csv_file):
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)

    # Append data to file
    with open(csv_file, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writerow({
            'command': data['command'],
            'user_time': data['user_time'],
            'system_time': data['system_time'],
            'elapsed_time': data['elapsed_time'],
            'cpu_percent': data['cpu_percent'],
            'max_resident_set_size': data['max_resident_set_size']
        })


for path in snakemake.input.time:
    append_to_csv(snakemake.output.time, parse_file(path))
    