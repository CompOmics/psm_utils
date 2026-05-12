"""psm_utils Streamlit-based web server."""

from _base import StreamlitPage


class StreamlitPageHome(StreamlitPage):
    """Streamlit page for the home section."""

    def _main_page(self):
        pass


if __name__ == "__main__":
    StreamlitPageHome()
