import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)
import io
B=OUTDIR
tpl=io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"report.template.html"),encoding="utf-8").read()
data=io.open(B+"/report_slim.json",encoding="utf-8").read()
assert "/*__DATA__*/" in tpl, "placeholder missing"
out=tpl.replace("/*__DATA__*/", data)
io.open(B+"/report.html","w",encoding="utf-8").write(out)
import os
print("report.html MB:", round(os.path.getsize(B+"/report.html")/1e6,2))
print("placeholder replaced:", "/*__DATA__*/" not in out)
print("has DATA designs:", '"designs"' in out, " all48:", '"all48"' in out)
