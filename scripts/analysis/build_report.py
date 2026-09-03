import io
B="/mnt/c/Users/thanh/AppData/Local/Temp/claude/C--Users-thanh-Documents/9a8e1c59-098d-4338-a1cc-e43c6d31654e/scratchpad"
tpl=io.open(B+"/report.template.html",encoding="utf-8").read()
data=io.open(B+"/report_slim.json",encoding="utf-8").read()
assert "/*__DATA__*/" in tpl, "placeholder missing"
out=tpl.replace("/*__DATA__*/", data)
io.open(B+"/report.html","w",encoding="utf-8").write(out)
import os
print("report.html MB:", round(os.path.getsize(B+"/report.html")/1e6,2))
print("placeholder replaced:", "/*__DATA__*/" not in out)
print("has DATA designs:", '"designs"' in out, " all48:", '"all48"' in out)
