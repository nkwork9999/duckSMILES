#!/usr/bin/env python3
"""Exercise documented numerical fixtures through the production exported C ABI."""
import argparse,ctypes,hashlib,json
from pathlib import Path
from check_jmp_continuation import equal,field,strict_json
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',required=True,type=Path);p.add_argument('--api',nargs='+',required=True);a=p.parse_args()
 lib=ctypes.CDLL(str(a.library.resolve()));run=lib.ds_analytics_run;run.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.POINTER(ctypes.c_double),ctypes.c_size_t,ctypes.c_size_t,ctypes.POINTER(ctypes.c_double),ctypes.c_size_t,ctypes.POINTER(ctypes.c_size_t)];run.restype=ctypes.c_void_p
 free=lib.ds_analytics_free;free.argtypes=[ctypes.c_void_p,ctypes.c_size_t];free.restype=None
 ops=[o for name in a.api for o in json.loads((ROOT/'tasks'/name).read_text())['operations']];count=0
 for op in ops:
  data=op['fixture']['data'];args=op['fixture']['args'];rows=len(data);columns=len(data[0]);assert rows and columns
  dv=(ctypes.c_double*(rows*columns))(*[v for row in data for v in row]);av=(ctypes.c_double*len(args))(*args);length=ctypes.c_size_t();name=op['operation'].encode()
  ptr=run(name,len(name),dv,rows,columns,av,len(args),ctypes.byref(length));assert ptr,(op['operation'],'C ABI returned null')
  try:result=strict_json(ctypes.string_at(ptr,length.value).decode())
  finally:free(ptr,length.value)
  for key,value in op['known_value'].items():equal(field(result,key),value,op['operation']+'.'+key,5e-5);count+=1
 print(json.dumps({'status':'PASS','layer':'production exported ds_analytics_run/free','operations':len(ops),'numeric_fields':count,'library_sha256':hashlib.sha256(a.library.read_bytes()).hexdigest()},indent=2))
if __name__=='__main__':main()
