## Changing port
...

## Visible accross your network
...

## SSH port forwarding from a remote machine
yes | pip install jupyter_http_over_ws
jupyter serverextension enable --py jupyter_http_over_ws;
...

--no-browser

## Cool theme
Grey what is this? Let's go 80s style!
```bash
jupyter labextension install @yeebc/jupyterlab_neon_theme
```

## Spelling iz Gud
```bash
yes | pip install jupyterlab-spellchecker
```

## Conda on Colab

```bash
!pip install -q condacolab
import condacolab
condacolab.install() #ignore message about session crashing, this is normal
```
