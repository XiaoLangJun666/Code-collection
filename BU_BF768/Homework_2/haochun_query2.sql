select DISTINCT g.gid, g.name
from targets t join gene g using (gid)
where t.mid in (SELECT mr2.mid FROM miRNA mr2  WHERE mr2.name like '%let-7b%') 
	and g.gid not in (select gid from targets t where t.mid in (select mr.mid from miRNA mr where mr.name like '%miR-18%'))


