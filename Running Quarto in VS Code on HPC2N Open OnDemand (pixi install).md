# Running Quarto in VS Code on HPC2N Open OnDemand (pixi install)

2026-11-09 · Mimmi Eriksson (support by Claude Opus 5.5)

## Background

Quarto installed with pixi (from conda-forge) only works inside an **activated** environment. Activation sets variables such as `QUARTO_DENO`, `QUARTO_SHARE_PATH` and `QUARTO_PANDOC` that tell the `quarto` script where its tools live. The VS Code Quarto extension calls `quarto` without activating anything, so the script looks in the wrong places and fails. The extension then reports "install Quarto". Every step below works around this.

In this doc, `<project>` stands for the absolute path to your project folder, the one with `pixi.toml`. Run `pwd` in that folder to get it, and replace `<project>` with it everywhere.

## Setup steps

Eight steps, in order. Each says what to do and why.

**1. Install Quarto and packages you need ****with pixi )**

Note: In this test Quarto was only tested using python, but one might want to have R too and more packages.

```bash
pixi add quarto pyyaml nbformat nbclient ipykernel
```

This gives a project-local Quarto plus its tools (deno, pandoc, sass, typst, esbuild) in `.pixi/envs/default/bin`. Quarto's Jupyter engine (`share/quarto/jupyter/jupyter.py`) imports PyYAML, nbformat and nbclient. Without them, rendering a `.qmd` with Python code fails with `ModuleNotFoundError`. If something is still missing, `pixi add jupyter` pulls in the full set.

**2. Create a wrapper script `<project>/bin/quarto`**

```bash
#!/usr/bin/env bash
ENV=<project>/.pixi/envs/default
export QUARTO_CONDA_PREFIX="$ENV"
export QUARTO_SHARE_PATH="$ENV/share/quarto"
export QUARTO_DENO="$ENV/bin/deno"
export QUARTO_DENO_DOM="$ENV/lib/deno_dom.so"
export DENO_DOM_PLUGIN="$ENV/lib/deno_dom.so"
export QUARTO_PANDOC="$ENV/bin/pandoc"
export QUARTO_DART_SASS="$ENV/bin/sass"
export QUARTO_ESBUILD="$ENV/bin/esbuild"
export QUARTO_TYPST="$ENV/bin/typst"
export QUARTO_TYPST_GATHER="$ENV/bin/typst-gather"
export PATH="$ENV/bin:$PATH"
exec "$ENV/bin/quarto" "$@"
```

The wrapper sets the variables that pixi's activation would set, then runs the real Quarto. It doesn't call pixi, so it behaves the same in a terminal and in the VS Code extension host, which has a very minimal environment. `QUARTO_DENO_DOM` and `DENO_DOM_PLUGIN` point to the correct file: the conda-forge package sets them to a path from its build machine (`/home/conda/feedstock_root/...`) that doesn't exist on HPC2N. Left alone, that breaks HTML post-processing.

The wrapper is created in `<project>/bin/quarto`, outside `.pixi`, only to keep a safe copy. Step 4 copies it into the environment, and `pixi add`, `pixi update` or `pixi install` can overwrite the copy there. The one in `<project>/bin/quarto` is what you restore from.

**3. Make the wrapper executable and test it**

```bash
chmod +x bin/quarto
bin/quarto check
env -i HOME="$HOME" PATH=/usr/bin:/bin bin/quarto --paths
```

`quarto check` should pass every step, including "Checking Jupyter engine render" with no traceback. A Dart Sass version note is only a warning. The `env -i` test runs the wrapper in an almost empty environment, like the extension host's. It should print two paths: the environment's `bin` and its `share/quarto` folder.

**4. Put the wrapper in place of the environment's own `quarto`**

```bash
B=<project>/.pixi/envs/default/bin
mv "$B/quarto" "$B/quarto.real"
cp <project>/bin/quarto "$B/quarto"
sed -i 's|exec "$ENV/bin/quarto"|exec "$ENV/bin/quarto.real"|' "$B/quarto"
chmod +x "$B/quarto"
env -i HOME="$HOME" PATH=/usr/bin:/bin "$B/quarto" --paths
```

`quarto --paths` reports `.pixi/envs/default/bin` as Quarto's bin folder, and the Quarto language server (LSP) runs the `quarto` in that folder directly, bypassing `quarto.path`. Without this step the LSP log shows `share/version: No such file`, followed by `ERR_MODULE_NOT_FOUND ... editor/tools/vs-code.mjs`, and the server stops after 5 crashes. Use absolute paths as above: the commands must not depend on the current folder.

`pixi add`, `pixi update` or `pixi install` can reinstall Quarto and restore the original script. If the LSP crash comes back after changing the environment, repeat this step.

**5. Point the Quarto extension at the wrapper**

In the **workspace** settings (`<project>/.vscode/settings.json`), set `quarto.path` to the wrapper's **absolute** path, `<project>/``.pixi/envs/default/``bin/quarto`. A relative path isn't resolved. Don't set it in User settings: other projects would then start the Quarto language server with this project's Quarto and log crashes. In projects where you don't use Quarto, disable the extension there (Extensions → Quarto → Disable (Workspace)).

**6. Turn off Python's automatic terminal activation**

```json
"python.terminal.activateEnvironment": false,
"quarto.usePipQuarto": false
```

With a pixi interpreter selected, the Python extension typed `pixi shell` into every new terminal. That starts a nested shell, which swallowed the `quarto preview` command the Quarto extension sent right after it. `quarto.usePipQuarto` off stops the Quarto extension from picking up the `quarto` that sits next to the selected Python.

**7. Set the Quarto variables for all VS Code terminals**

In the workspace settings (`.vscode/settings.json`):

```json
"terminal.integrated.env.linux": {
  "QUARTO_SHARE_PATH": "<project>/.pixi/envs/default/share/quarto",
  "QUARTO_DENO": "<project>/.pixi/envs/default/bin/deno",
  "QUARTO_DENO_DOM": "<project>/.pixi/envs/default/lib/deno_dom.so",
  "QUARTO_PANDOC": "<project>/.pixi/envs/default/bin/pandoc",
  "QUARTO_DART_SASS": "<project>/.pixi/envs/default/bin/sass",
  "QUARTO_ESBUILD": "<project>/.pixi/envs/default/bin/esbuild",
  "QUARTO_TYPST": "<project>/.pixi/envs/default/bin/typst",
  "QUARTO_TYPST_GATHER": "<project>/.pixi/envs/default/bin/typst-gather",
  "QUARTO_CONDA_PREFIX": "<project>/.pixi/envs/default"
}
```

The Preview terminal ignored the wrapper. An extension (probably Python Environments) puts `.pixi/envs/default/bin` on PATH, so plain `quarto` resolved to the environment's copy. With these variables set in every VS Code terminal, that copy works however it's called. Python and Jupyter from pixi come in through the same PATH entry. The values are what `pixi run env | grep -i quarto` prints, except `QUARTO_DENO_DOM`, corrected as in step 2.

**8. Reload and use Preview**

Run "Developer: Reload Window", close all old terminals (they keep their old environment), then click Preview. The rendered `.html` appears next to the `.qmd`.

## Known limitation: Quarto's own preview doesn't load

Rendering works, but the preview page from Quarto's Preview button doesn't load, so no side panel opens. Quarto prints a link through code-server's port proxy, e.g.:

```
Browse at https://portal.hpc2n.umu.se/rnode/<node>/<port>/proxy/<preview-port>/<file>.html
Listening on http://127.0.0.1:<preview-port>/
```

`<preview-port>` is the port Quarto picks for its preview server, and it changes on each run. `<node>`, `<port>` and `<file>` are explained under "Getting a preview side panel" below.

That link gives a white page with "Not found". What was tested:

| Test | Result | Meaning |
| --- | --- | --- |
| `curl -sI http://127.0.0.1:<preview-port>/<file>.html` on the node | `200 OK` | Quarto's preview server works locally |
| code-server flags (`ps -o args= -u $USER \| grep [c]ode-server`) | No `--disable-proxy` | code-server's port proxy is on |
| `python -m http.server 8765`, opened at `.../<port>/proxy/8765/<file>.html` | Page displays | OnDemand → code-server → port works |
| Quarto's link via `.../absproxy/<preview-port>/...` | Not found | Expected: Quarto needs the prefix stripped (`proxy`) |

So Quarto's preview server is the part that rejects proxied requests, probably because of a header the proxy adds. The cause isn't known yet. The preview side panel options below work around it.

## Getting a preview side panel

Both options serve the rendered `.html` with a plain web server, which works through the proxy, and show it in a VS Code panel. Render with the Render button or `quarto render <file>.qmd`.

**Option 1: Simple Browser (built in, no install)**

1. In a terminal in `<project>`, start a server and leave it running: `python -m http.server 8765`
2. Open the Command Palette (Ctrl+Shift+P), run **"Simple Browser: Show"**, and paste `https://portal.hpc2n.umu.se/rnode/<node>/<port>/proxy/8765/<file>.html`
3. Drag the panel next to your code. Click its reload button after each render.

Fill in the placeholders for your session and file:

- `<node>`: the compute node your OnDemand session runs on, e.g. `b-cn1701.hpc2n.umu.se`.
- `<port>`: code-server's port, e.g. `44451`.
- `<file>`: the path to the rendered file, relative to the folder where you started the server, e.g. `notebooks/test`.

`<node>` and `<port>` change with each OnDemand session. Copy them from the browser address bar, which shows `.../rnode/<node>/<port>/...` while VS Code is open.

**Option 2: Live Preview extension (auto-reload)**

1. Install Microsoft's **Live Preview** extension (`ms-vscode.live-server`) from the Extensions panel. Check that it is available in code-server's marketplace (Open VSX).
2. Right-click `<file>.html` and choose **"Show Preview"**.
3. Render the `.qmd`. The panel reloads by itself when the `.html` changes.

Live Preview runs its own server and goes through code-server's port proxy, which is shown to work above. Option 2 hasn't been tested on HPC2N yet.

## Troubleshooting

Match the error to its cause and fix:

| Symptom | Cause | Fix |
| --- | --- | --- |
| "Install Quarto" when clicking Preview | Extension can't run Quarto | Steps 2–5 |
| `cat: .../share/version: No such file or directory` | `quarto` run without activation | Steps 2, 4 and 7 |
| `.../bin/tools/x86_64/deno: No such file or directory` | `quarto` run without activation | Steps 2, 4 and 7 |
| Terminal shows `pixi shell ...` then nothing happens | Python extension auto-activation | Step 6 |
| `ModuleNotFoundError: No module named 'yaml'` (or `nbformat`) | Jupyter packages missing | Step 1 |
| Quarto LSP server crashed 5 times, with `ERR_MODULE_NOT_FOUND ... editor/tools/vs-code.mjs` | LSP runs the environment's `quarto` directly | Step 4, then reload |
| Same LSP crash in a project that doesn't use Quarto | `quarto.path` set in User settings | Step 5 |
| LSP crash returns after `pixi add` / `pixi update` | Original `quarto` script restored | Repeat step 4 |
| Preview link shows "Not found", no side panel | Quarto's preview server rejects proxied requests | Use a side panel option above |

Useful checks:

- `type -a quarto` in a new VS Code terminal shows which `quarto` runs.
- `pixi run env | grep -i quarto` shows the values activation sets.
- `bin/quarto check` tests the whole setup outside VS Code.
- `env -i HOME="$HOME" PATH=/usr/bin:/bin <path to quarto> --paths` tests a `quarto` the way the extension host runs it.
- The "Quarto" channel in the Output tab shows why the language server crashed. Read the lines above the crash message.

Some things turned out not to matter. A first wrapper using `pixi run` worked by hand but was never used by Preview, because of the PATH issue in step 7. A wrapper using `pixi shell-hook` also worked, but didn't fix the language server, because of the bin path in step 4. The different filesystems (`/home` vs `/pfs`) didn't break anything: `bin/quarto` is a normal file in this build, not a symlink.
