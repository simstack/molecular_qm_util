import asyncio
import logging
import warnings
from pathlib import Path

import pytest
import pytest_asyncio

from simstack.core.context import context
from simstack.tables.model_table import make_model_table
from simstack.tables.node_table import make_node_table
from simstack.util.project_root_finder import find_project_root

# Suppress pymongo logs/warnings
logging.getLogger("pymongo").setLevel(logging.WARNING)
# Suppress motor logs
logging.getLogger("motor").setLevel(logging.WARNING)

# Suppress Pydantic deprecation warnings
warnings.filterwarnings("ignore", message=".*json_encoders.*")
warnings.filterwarnings("ignore", message=".*model_computed_fields.*")
warnings.filterwarnings("ignore", message=".*model_fields.*")
try:
    from pymongo.errors import PyMongoDeprecationWarning
    warnings.filterwarnings("ignore", category=PyMongoDeprecationWarning)
except ImportError:
    pass


def pytest_addoption(parser):
    parser.addoption("--gather", action="store_true", default=False, help="gather result data for tests")


@pytest.fixture
def gather(request):
    return request.config.getoption("--gather")


@pytest.fixture(scope="session")
def event_loop():
    """Create an instance of the default event loop for the test session."""
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


@pytest_asyncio.fixture(scope="session", autouse=True, loop_scope="session")
async def initialized_context():
    """
    Session-based async fixture that initializes the simstack context.
    """
    project_root = find_project_root()

    await context.initialize()

    dirs = [Path(__file__).resolve().parents[1] / "simstack-confgen"]

    await make_model_table(
        context.db,
        dirs=dirs,
        drops="src",
        clear=False,
        project_root=project_root,
        ignore_entrypoints=False,
    )
    await make_node_table(
        context.db,
        dirs=dirs,
        drops="src",
        clear=False,
        project_root=project_root,
        ignore_entrypoints=False,
    )

    await context.refresh_mappings()

    yield context

    # Cleanup
    if context.initialized:
        if hasattr(context, "_db") and context._db:
            await context._db.close()
            context._db = None

        context._initialized = False
        context._model_mappings = None
        context._node_mappings = None
        context._resource_config = None
        context._config = None
