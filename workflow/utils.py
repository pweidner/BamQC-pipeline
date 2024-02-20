import subprocess

def extract_alignment_metrics(tsv_file):
    # Use zgrep and cut to extract alignment summary metrics
    grep_process = subprocess.Popen(['zgrep', '^ME', tsv_file], stdout=subprocess.PIPE)
    cut_process = subprocess.Popen(['cut', '-f', '2-'], stdin=grep_process.stdout, stdout=subprocess.PIPE)
    grep_process.stdout.close()  # Allow grep process to receive a SIGPIPE if cut exits
    alignment_metrics = {}
    with cut_process.stdout as f:
        for line_bytes in f:
            line = line_bytes.decode('utf-8').strip()  # Decode bytes to string
            fields = line.split("\t")
            sample = fields[0]
            metrics = fields[1:]
            alignment_metrics[sample] = metrics
    return alignment_metrics
