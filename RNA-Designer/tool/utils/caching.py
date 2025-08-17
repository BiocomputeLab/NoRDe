import pickle
from pathlib import Path
import atexit

class SequenceCache:
    def __init__(self, cache_file):
        self.cache_file = Path(cache_file)
        self.cache = self._load_cache()
        
    def _load_cache(self):
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'rb') as f:
                    return pickle.load(f)
            except (pickle.PickleError, EOFError):
                return {}
        return {}
    
    def save(self):
        with open(self.cache_file, 'wb') as f:
            pickle.dump(self.cache, f)
            
    def get(self, key):
        return self.cache.get(key)
    
    def set(self, key, value):
        self.cache[key] = value

class CacheManager:
    def __init__(self):
        self.rnafold_cache = SequenceCache("output/rnafold_cache.pkl")
        self.gc_content_cache = SequenceCache("output/gc_content_cache.pkl")
        atexit.register(self.rnafold_cache.save)
        atexit.register(self.gc_content_cache.save)