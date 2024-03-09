import subprocess

program_path = "../poisson.exe"


def average_time(n, *args):
    results_sum = 0
    for i in range(n):
        result = subprocess.run([program_path, *args, "--time-only"], capture_output=True)
        if result.returncode:
            print("Error in program")
            print(result.stderr)
        results_sum += float(result.stdout)
    return results_sum / n


def preset_1(*args):
    print(average_time(10, "-n 98", "-c 20", "-t \"task_2\"", *args))


preset_1("-a 1")
preset_1("-a 3", "-j 1")
preset_1("-a 3", "-j 2")
preset_1("-a 3", "-j 3")
preset_1("-a 3", "-j 4")
preset_1("-a 3", "-j 5")
preset_1("-a 3", "-j 6")
preset_1("-a 3", "-j 12")
preset_1("-a 6", "-j 1")
preset_1("-a 6", "-j 2")
preset_1("-a 6", "-j 3")
preset_1("-a 6", "-j 4")
preset_1("-a 6", "-j 5")
preset_1("-a 6", "-j 6")
preset_1("-a 6", "-j 12")


def preset_2(*args):
    print(average_time(3, "-n 498", "-c 50", "-t \"task_2\"", *args))


preset_2("-a 1")
preset_2("-a 3", "-j 1")
preset_2("-a 3", "-j 2")
preset_2("-a 3", "-j 3")
preset_2("-a 3", "-j 4")
preset_2("-a 3", "-j 5")
preset_2("-a 3", "-j 6")
preset_2("-a 3", "-j 12")
preset_2("-a 6", "-j 1")
preset_2("-a 6", "-j 2")
preset_2("-a 6", "-j 3")
preset_2("-a 6", "-j 4")
preset_2("-a 6", "-j 5")
preset_2("-a 6", "-j 6")
preset_2("-a 6", "-j 12")


def preset_3(*args):
    print(average_time(3, "-n 498", "-c 100", "-t \"task_2\"", *args))


preset_3("-a 6", "-j 1")
preset_3("-a 6", "-j 2")
preset_3("-a 6", "-j 3")
preset_3("-a 6", "-j 4")
preset_3("-a 6", "-j 5")
preset_3("-a 6", "-j 6")
preset_3("-a 6", "-j 12")


def preset_4(*args):
    print(average_time(1, "-n 998", "-c 50", "-e 0.01", "-t \"task_2\"", *args))


preset_4("-a 1")
preset_4("-a 3", "-j 1")
preset_4("-a 3", "-j 2")
preset_4("-a 3", "-j 3")
preset_4("-a 3", "-j 5")
preset_4("-a 3", "-j 4")
preset_4("-a 3", "-j 6")
preset_4("-a 3", "-j 12")
preset_4("-a 6", "-j 1")
preset_4("-a 6", "-j 2")
preset_4("-a 6", "-j 3")
preset_4("-a 6", "-j 4")
preset_4("-a 6", "-j 5")
preset_4("-a 6", "-j 6")
preset_4("-a 6", "-j 12")


def preset_5(*args):
    print(average_time(1, "-n 3998", "-c 100", "-e 0.1", "-t \"task_2\"", *args))


preset_5("-a 1")
preset_5("-a 3", "-j 1")
preset_5("-a 3", "-j 2")
preset_5("-a 3", "-j 3")
preset_5("-a 3", "-j 5")
preset_5("-a 3", "-j 4")
preset_5("-a 3", "-j 6")
preset_5("-a 3", "-j 12")
preset_5("-a 6", "-j 1")
preset_5("-a 6", "-j 2")
preset_5("-a 6", "-j 3")
preset_5("-a 6", "-j 4")
preset_5("-a 6", "-j 5")
preset_5("-a 6", "-j 6")
preset_5("-a 6", "-j 12")
