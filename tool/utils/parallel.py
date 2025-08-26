import logging
import multiprocessing
from functools import partial
from typing import Callable, Iterable, Any, List
from tool.config import Config

logger = logging.getLogger(__name__)

class ParallelProcessor:
    """
    Handles parallel processing of tasks with configurable number of cores.
    
    Provides both map and starmap functionality with automatic core detection.
    """
    
    @staticmethod
    def get_available_cores() -> int:
        """Get number of available CPU cores with safety margin."""
        return max(1, multiprocessing.cpu_count() - 1)
    
    @classmethod
    def parallel_map(cls, 
                    func: Callable, 
                    iterable: Iterable, 
                    n_cores: int = None,
                    chunksize: int = None) -> List[Any]:
        """
        Parallel map operation.
        
        Args:
            func: Function to apply
            iterable: Items to process
            n_cores: Number of cores to use (default: auto-detect)
            chunksize: Chunk size for workload distribution
            
        Returns:
            List of results
        """
        n_cores = n_cores or Config.NUM_CORES or cls.get_available_cores()
        
        try:
            with multiprocessing.Pool(processes=n_cores) as pool:
                return pool.map(func, iterable, chunksize=chunksize)
        except Exception as e:
            logger.error(f"Parallel map error: {e}")
            logger.info("Falling back to sequential processing")
            return [func(item) for item in iterable]
    
    @classmethod
    def parallel_starmap(cls,
                        func: Callable,
                        iterable: Iterable[tuple],
                        n_cores: int = None,
                        chunksize: int = None) -> List[Any]:
        """
        Parallel starmap operation (for functions with multiple arguments).
        
        Args:
            func: Function to apply
            iterable: Tuples of arguments to process
            n_cores: Number of cores to use (default: auto-detect)
            chunksize: Chunk size for workload distribution
            
        Returns:
            List of results
        """
        n_cores = n_cores or Config.NUM_CORES or cls.get_available_cores()
        
        try:
            with multiprocessing.Pool(processes=n_cores) as pool:
                return pool.starmap(func, iterable, chunksize=chunksize)
        except Exception as e:
            logger.error(f"Parallel starmap error: {e}")
            logger.info("Falling back to sequential processing")
            return [func(*args) for args in iterable]
    
    @classmethod
    def parallel_apply(cls,
                      func: Callable,
                      constant_args: dict,
                      variable_args: Iterable[dict],
                      n_cores: int = None) -> List[Any]:
        """
        Parallel application of function with mixed constant and variable arguments.
        
        Args:
            func: Function to apply
            constant_args: Dictionary of arguments that stay the same for all calls
            variable_args: Iterable of dictionaries with arguments that vary
            n_cores: Number of cores to use (default: auto-detect)
            
        Returns:
            List of results
        """
        n_cores = n_cores or Config.NUM_CORES or cls.get_available_cores()
        partial_func = partial(cls._apply_args, func, constant_args)
        
        try:
            with multiprocessing.Pool(processes=n_cores) as pool:
                return pool.map(partial_func, variable_args)
        except Exception as e:
            logger.error(f"Parallel apply error: {e}")
            logger.info("Falling back to sequential processing")
            return [partial_func(args) for args in variable_args]
    
    @staticmethod
    def _apply_args(func: Callable, constant_args: dict, variable_args: dict) -> Any:
        """Helper method to merge constant and variable arguments."""
        kwargs = {**constant_args, **variable_args}
        return func(**kwargs)