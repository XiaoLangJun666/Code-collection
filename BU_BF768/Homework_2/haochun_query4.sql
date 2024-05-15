select p1.gid as gene_id_1,p1.name as gene_name_1,p2.gid as gene_id_2,p2.name as gene_name_2,count(p1.mid) as count_pairs
from (select g1.gid,g1.name,t1.mid from targets t1 join gene g1 using (gid)) as p1
join (select g2.gid,g2.name,t2.mid from targets t2 join gene g2 using (gid)) as p2 on p1.mid = p2.mid and p1.gid < p2.gid
group by p1.gid, p2.gid
ORDER BY count_pairs DESC
limit 5



