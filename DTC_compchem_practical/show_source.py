def display_source(function) -> None:
    """
    Display the source code of a function or class.

    Notebook only. Won't work for C code.
    """
    import inspect
    from IPython.display import HTML, display
    from pygments import highlight
    from pygments.lexers import PythonLexer
    from pygments.formatters import HtmlFormatter
    code:str = inspect.getsource(function)
    html:str = highlight(code, PythonLexer(), HtmlFormatter(style='colorful'))
    stylesheet:str = f"<style>{HtmlFormatter().get_style_defs('.highlight')}</style>"
    display(HTML(f"{stylesheet}{html}"))