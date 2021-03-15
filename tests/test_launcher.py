import os
import tempfile
import dask.distributed
import swot_simulator.launcher
import swot_simulator.settings

ROOT = os.path.dirname(os.path.abspath(__file__))


def test_launcher():
    with tempfile.TemporaryDirectory() as tmpdir:
        parameters = swot_simulator.settings.Parameters.load_default()
        parameters.ephemeris = os.path.join(
            ROOT, "..", "data", "ephemeris_calval_june2015_ell.txt")
        parameters.nadir = True
        parameters.working_directory = tmpdir
        cluster = dask.distributed.LocalCluster()
        client = dask.distributed.Client(cluster)
        client.wait_for_workers(1)
        swot_simulator.launcher.launch(client, parameters, None)
        client.close()
        cluster.close()
