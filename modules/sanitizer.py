import html
import re
from typing import Any, Union, Dict

class OutputSanitizer:
    @staticmethod
    def sanitize_html(text: str) -> str:
        """Sanitize HTML to prevent XSS."""
        return html.escape(text)
    
    @staticmethod
    def sanitize_filename(filename: str) -> str:
        """Sanitize filenames to prevent path traversal and command injection."""
        # Remove any non-alphanumeric characters except for dots and underscores
        sanitized = re.sub(r'[^a-zA-Z0-9._-]', '', filename)
        # Ensure no leading periods (hidden files) or trailing spaces
        sanitized = sanitized.lstrip('.').strip()
        return sanitized
    
    @staticmethod
    def sanitize_error_message(message: str) -> str:
        """Sanitize error messages to prevent information disclosure."""
        # Remove potential system paths
        message = re.sub(r'(?i)([A-Z]:\\|/var/|/etc/|/usr/|/home/|~/).*?(?=\s|$)', '[PATH]', message)
        # Remove potential stack traces
        message = re.sub(r'File.*?line \d+.*?\n', '[TRACE]\n', message)
        # Remove potential SQL queries
        message = re.sub(r'(SELECT|INSERT|UPDATE|DELETE|DROP).*?;', '[SQL]', message, flags=re.I)
        return message
    
    @staticmethod
    def sanitize_json_output(data: Any) -> Any:
        """Recursively sanitize JSON output."""
        if isinstance(data, dict):
            return {
                OutputSanitizer.sanitize_html(str(k)): 
                OutputSanitizer.sanitize_json_output(v)
                for k, v in data.items()
            }
        elif isinstance(data, list):
            return [OutputSanitizer.sanitize_json_output(item) for item in data]
        elif isinstance(data, str):
            return OutputSanitizer.sanitize_html(data)
        else:
            return data
    
    @staticmethod
    def sanitize_config_output(config: Dict[str, Any]) -> Dict[str, Any]:
        """Sanitize configuration output while preserving structure."""
        sanitized = {}
        for key, value in config.items():
            if isinstance(value, dict):
                sanitized[OutputSanitizer.sanitize_html(key)] = OutputSanitizer.sanitize_config_output(value)
            elif isinstance(value, list):
                sanitized[OutputSanitizer.sanitize_html(key)] = [
                    OutputSanitizer.sanitize_html(str(item)) if isinstance(item, str) else item
                    for item in value
                ]
            elif isinstance(value, str):
                sanitized[OutputSanitizer.sanitize_html(key)] = OutputSanitizer.sanitize_html(value)
            else:
                sanitized[OutputSanitizer.sanitize_html(key)] = value
        return sanitized