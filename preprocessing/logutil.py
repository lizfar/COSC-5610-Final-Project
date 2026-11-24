"""Colored progress messages for every preprocessing step."""

from termcolor import cprint


def step(message):
    cprint(f"[STEP] {message}", "cyan")


def ok(message):
    cprint(f"[OK] {message}", "green")


def warn(message):
    cprint(f"[WARN] {message}", "yellow")


def fail(message):
    cprint(f"[ERROR] {message}", "red")


class StepLog:
    """Collect the same sentences that are printed, for the written report."""

    def __init__(self, title):
        self.title = title
        self.lines = []
        step(title)
        self.lines.append(title)

    def add(self, message):
        step(message)
        self.lines.append(message)

    def add_warn(self, message):
        warn(message)
        self.lines.append(f"WARNING: {message}")

    def add_ok(self, message):
        ok(message)
        self.lines.append(message)
