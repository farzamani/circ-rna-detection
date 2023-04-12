import os
import csv
import datetime


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
            result['elapsed_time'] = line.split(": ")[1].strip()
        elif line.startswith("\tMaximum resident set size"):
            result['max_resident_set_size'] = int(line.split(": ")[1])
    
    return result


def compile_files(input_files, output_file):
    max_cpu_percent = 0
    max_resident_set_size = 0
    total_user_time = 0
    total_system_time = 0
    total_elapsed_time = datetime.timedelta()

    for path in input_files:
        result = parse_file(path)
        command = result.get('command', 0)
        max_cpu_percent = max(max_cpu_percent, result.get('cpu_percent', 0))
        max_resident_set_size = max(max_resident_set_size, result.get('max_resident_set_size', 0))
        total_user_time += result.get('user_time', 0)
        total_system_time += result.get('system_time', 0)
        elapsed_time_str = result.get('elapsed_time', '0:00')
        elapsed_time_parts = elapsed_time_str.split(':')
        if len(elapsed_time_parts) == 2:
            elapsed_time_parts = ['0'] + elapsed_time_parts
        elapsed_time_seconds = float(elapsed_time_parts[2])
        elapsed_time_seconds = int(elapsed_time_seconds)
        elapsed_time = datetime.timedelta(hours=int(elapsed_time_parts[0]),
                                        minutes=int(elapsed_time_parts[1]),
                                        seconds=elapsed_time_seconds)
        total_elapsed_time += elapsed_time

    with open(output_file, 'w') as f:
        f.write(f"\tCommand being timed: {command}\n")
        f.write(f"\tUser time (seconds): {total_user_time:.2f}\n")
        f.write(f"\tSystem time (seconds): {total_system_time:.2f}\n")
        elapsed_time_str = str(total_elapsed_time).split('.')[0]
        f.write(f"\tElapsed (wall clock) time (h:mm:ss or m:ss): {elapsed_time_str}\n")
        f.write(f"\tPercent of CPU this job got: {max_cpu_percent}%\n")
        f.write(f"\tMaximum resident set size (kbytes): {max_resident_set_size}\n")


compile_files(snakemake.input, snakemake.output.time)
