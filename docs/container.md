# Container usage

Build locally from the repository root:

```bash
docker build -t proteovault:local .
```

Run the CLI with the directories containing both the pipeline output and database
destination mounted into the container:

```bash
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v /path/to/pipeline-output:/input:ro \
    -v /path/to/results:/output \
    proteovault:local \
    proteovault build --pipeline-output /input --output /output/proteovault.duckdb

docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v /path/to/results:/output:ro \
    proteovault:local \
    proteovault validate /output/proteovault.duckdb
```

The image has no long-running default process. Supply a ProteoVault command to
`docker run`, or supply an explicit long-running command if using `docker exec`.
