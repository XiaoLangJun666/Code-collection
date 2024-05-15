select g.gid, g.name, count(t.mid) as miRNA_count
from targets t join gene g using (gid)
where g.name not like 'A%'
GROUP BY g.gid
ORDER BY miRNA_count DESC
limit 10


