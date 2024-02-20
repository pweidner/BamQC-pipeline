import subprocess

def extract_metrics(tsv_file, row_id):
    # Use zgrep and cut to extract metrics based on row_id
    grep_process = subprocess.Popen(['zgrep', f'^{row_id}', tsv_file], stdout=subprocess.PIPE)
    cut_process = subprocess.Popen(['cut', '-f', '2-'], stdin=grep_process.stdout, stdout=subprocess.PIPE)
    grep_process.stdout.close()  # Allow grep process to receive a SIGPIPE if cut exits
    metrics = {}
    with cut_process.stdout as f:
        for line_bytes in f:
            line = line_bytes.decode('utf-8').strip()  # Decode bytes to string
            fields = line.split("\t")
            sample = fields[0]
            data = fields[1:]
            metrics[sample] = data
    return metrics
