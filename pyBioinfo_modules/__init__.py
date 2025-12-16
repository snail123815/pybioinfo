import logging

# Configure logging once when the package is first imported
if not logging.getLogger().hasHandlers():
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
