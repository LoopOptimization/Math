#!/usr/bin/env python3                                                                   

import json
from pathlib import Path

files = [
  "./build/bench/compile_commands.json",
  "./build/no-san/compile_commands.json",
]

commands = {}
for f in files:
    for entry in json.loads(Path(f).read_text()):
        commands[entry["file"]] = entry

Path("./compile_commands.json").write_text(json.dumps(list(commands.values()), indent=2))  

