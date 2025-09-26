import os
import re
from pathlib import Path
from typing import Optional, List
import magic  # python-magic for secure file type detection
from dataclasses import dataclass

@dataclass
class SecurityConfig:
    ALLOWED_EXTENSIONS = {'json'}
    MAX_FILE_SIZE = 1024 * 1024  # 1MB
    MAX_JSON_DEPTH = 20
    ALLOWED_JSON_TYPES = {str, int, float, bool, list, dict, None}
    PATH_TRAVERSAL_PATTERN = re.compile(r'\.{2}|/|\\')
    
class SecurityValidator:
    @staticmethod
    def is_safe_json_content(content: dict, depth: int = 0) -> bool:
        """Validate JSON content for security issues."""
        if depth > SecurityConfig.MAX_JSON_DEPTH:
            raise ValueError("JSON nesting too deep")
            
        if isinstance(content, dict):
            return all(
                isinstance(k, str) and 
                SecurityValidator.is_safe_json_content(v, depth + 1)
                for k, v in content.items()
            )
        elif isinstance(content, list):
            return all(
                SecurityValidator.is_safe_json_content(item, depth + 1)
                for item in content
            )
        else:
            return type(content) in SecurityConfig.ALLOWED_JSON_TYPES
            
    @staticmethod
    def is_safe_filename(filename: str) -> bool:
        """Check if filename is safe."""
        return '.' in filename and \
            filename.rsplit('.', 1)[1].lower() in SecurityConfig.ALLOWED_EXTENSIONS and \
            not SecurityConfig.PATH_TRAVERSAL_PATTERN.search(filename)
            
    @staticmethod
    def get_secure_file_type(file_path: str) -> Optional[str]:
        """Get secure file type using python-magic."""
        try:
            mime = magic.Magic(mime=True)
            file_type = mime.from_file(file_path)
            return file_type if file_type == 'application/json' else None
        except Exception:
            return None
            
    @staticmethod
    def sanitize_ethnicity_name(name: str) -> str:
        """Sanitize ethnicity name for safe file operations."""
        # Remove any characters that aren't alphanumeric, space, or hyphen
        sanitized = re.sub(r'[^a-zA-Z0-9\s-]', '', name)
        # Convert spaces to underscores and ensure no double hyphens
        sanitized = re.sub(r'[\s-]+', '_', sanitized)
        return sanitized.lower()
        
    @staticmethod
    def is_safe_file_size(file_path: str) -> bool:
        """Check if file size is within allowed limits."""
        return os.path.getsize(file_path) <= SecurityConfig.MAX_FILE_SIZE
        
    @staticmethod
    def secure_path_join(*paths: List[str]) -> str:
        """Securely join paths preventing directory traversal."""
        base = Path(paths[0]).resolve()
        try:
            joined = base.joinpath(*paths[1:]).resolve()
            if not str(joined).startswith(str(base)):
                raise ValueError("Path traversal detected")
            return str(joined)
        except (TypeError, ValueError):
            raise ValueError("Invalid path")