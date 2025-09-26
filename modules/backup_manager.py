import shutil
import os
from pathlib import Path
from datetime import datetime
import json
from typing import Optional, List
import hashlib
import gzip

class BackupManager:
    def __init__(self, config_dir: str, backup_dir: str = "backups"):
        self.config_dir = Path(config_dir)
        self.backup_dir = Path(backup_dir)
        self.backup_dir.mkdir(exist_ok=True)
        
    def calculate_file_hash(self, file_path: Path) -> str:
        """Calculate SHA-256 hash of a file."""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()
    
    def create_backup(self) -> str:
        """Create a compressed backup of all configurations."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = self.backup_dir / f"config_backup_{timestamp}.tar.gz"
        
        # Create manifest
        manifest = {
            "timestamp": timestamp,
            "files": {}
        }
        
        # Create temporary directory for backup
        temp_dir = self.backup_dir / f"temp_{timestamp}"
        temp_dir.mkdir()
        
        try:
            # Copy files and calculate hashes
            for config_file in self.config_dir.glob("*.json"):
                file_hash = self.calculate_file_hash(config_file)
                manifest["files"][config_file.name] = {
                    "hash": file_hash,
                    "size": config_file.stat().st_size
                }
                shutil.copy2(config_file, temp_dir)
            
            # Save manifest
            manifest_path = temp_dir / "manifest.json"
            with open(manifest_path, 'w') as f:
                json.dump(manifest, f, indent=2)
            
            # Create compressed archive
            shutil.make_archive(
                str(backup_path).rsplit('.', 1)[0],
                'gztar',
                temp_dir
            )
            
        finally:
            # Clean up temporary directory
            shutil.rmtree(temp_dir)
        
        return str(backup_path)
    
    def restore_backup(self, backup_path: str, validate: bool = True) -> bool:
        """Restore configurations from a backup."""
        backup_path = Path(backup_path)
        if not backup_path.exists():
            raise FileNotFoundError("Backup file not found")
        
        # Create temporary directory for restoration
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        temp_dir = self.backup_dir / f"restore_{timestamp}"
        temp_dir.mkdir()
        
        try:
            # Extract archive
            shutil.unpack_archive(backup_path, temp_dir)
            
            # Load and verify manifest
            manifest_path = temp_dir / "manifest.json"
            with open(manifest_path) as f:
                manifest = json.load(f)
            
            if validate:
                # Verify file hashes
                for filename, file_info in manifest["files"].items():
                    file_path = temp_dir / filename
                    if not file_path.exists():
                        raise ValueError(f"Missing file in backup: {filename}")
                    
                    current_hash = self.calculate_file_hash(file_path)
                    if current_hash != file_info["hash"]:
                        raise ValueError(f"Hash mismatch for file: {filename}")
            
            # Backup current configurations before restore
            self.create_backup()
            
            # Copy files to config directory
            for file_path in temp_dir.glob("*.json"):
                if file_path.name != "manifest.json":
                    shutil.copy2(file_path, self.config_dir)
            
            return True
            
        finally:
            # Clean up temporary directory
            shutil.rmtree(temp_dir)
            
    def list_backups(self) -> List[dict]:
        """List available backups with their information."""
        backups = []
        for backup_file in self.backup_dir.glob("*.tar.gz"):
            backup_info = {
                "filename": backup_file.name,
                "size": backup_file.stat().st_size,
                "created": datetime.fromtimestamp(backup_file.stat().st_mtime)
                            .strftime("%Y-%m-%d %H:%M:%S")
            }
            backups.append(backup_info)
        return sorted(backups, key=lambda x: x["created"], reverse=True)
    
    def cleanup_old_backups(self, keep_count: int = 10) -> int:
        """Remove old backups keeping only the specified number of recent ones."""
        backups = self.list_backups()
        if len(backups) <= keep_count:
            return 0
        
        removed_count = 0
        for backup in backups[keep_count:]:
            backup_path = self.backup_dir / backup["filename"]
            backup_path.unlink()
            removed_count += 1
        
        return removed_count