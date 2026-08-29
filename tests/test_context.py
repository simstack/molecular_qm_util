import pytest
from simstack.core.context import context


@pytest.mark.asyncio
async def test_initialized_context():
    assert context.initialized is True
    assert context.db is not None
    assert context.model_mappings is not None
    assert context.node_mappings is not None
