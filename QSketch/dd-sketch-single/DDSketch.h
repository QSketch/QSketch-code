#ifndef _DD_SKETCH_H_
#define _DD_SKETCH_H_

#include <bits/stdc++.h>
#include <hash.h>

#define inf 2147483647
#define log_base 1.1

class Slot {
public:
	Slot() {}
	Slot(uint32_t s): start(s), end(s + 1), f(1) {}
	~Slot() {}
	void create(Slot &t) {
		start = t.start;
		end = t.end;
		f = 1;
	}
	bool insert(uint32_t T) {
		if (start <= T && T < end) {
			f++;
			return true;
		}
		return false;
	} 
	void merge(Slot &t) {
		assert(t.start >= end);
		end = t.end;
		f += t.f;
	}
	bool operator < (Slot &t) {
		//printf("%d %d %d %d\n", start, end, t.start, t.end);
		assert((end <= t.start || start >= t.end));
		return start < t.start;
	}
	uint32_t start;
	uint32_t end;
	uint32_t f;
};

class Bucket {
public:
	Bucket() {}
	~Bucket() {
		delete[] slot;
	}
	void init(uint32_t c) {
		cell_num = c;
		slot = new Slot [cell_num];
		memset(slot, 0, cell_num * sizeof(Slot));
	}
	void insert (uint32_t T) {
		f++;
		for (int i = 0; i < slot_num; ++i) {
			if (slot[i].insert(T)) {
				std::sort(slot, slot + slot_num);
				return;
			}
		}
		// fail to insert, create a new slot
		/*for (int i = 0; i < slot_num; i++)
		{
			printf("%d %d    ", slot[i].start, slot[i].end);
		}
		printf("%d", T);
		printf("\n");*/
		Slot s = Slot(T);
		if (slot_num < cell_num) {
			slot[slot_num].create(s);
			slot_num++;
			std::sort(slot, slot + slot_num);
		}
		else {
			uint32_t mn = inf, mn1, mn2;
			int min_i = 0, p = cell_num;
			for (int i = 0; i < cell_num; i++)
			{
				if (s < slot[i]) 
				{
					p = i;
					break;
				}
			}
			if (p > 0) mn1 = slot[p - 1].f + s.f;
			else mn1 = inf;
			if (p < cell_num) mn2 = slot[p].f + s.f;
			else mn2 = inf;
			for (int i = 0; i < cell_num - 1; i++)
			{
				if (i == p - 1) continue;
				if (mn > slot[i].f + slot[i + 1].f)
				{
					mn = slot[i].f + slot[i + 1].f;
					min_i = i;
				}
			} 
			if (mn <= mn1 && mn <= mn2)
			{
				slot[min_i].merge(slot[min_i + 1]);
				slot[min_i + 1] = s;
				std::sort(slot, slot + cell_num);
			}
			else if (mn1 <= mn2)
			{
				slot[p - 1].merge(s);
			}
			else
			{
				s.merge(slot[p]);
				slot[p] = s;
			}
			/*// merge
			if (s < slot[0]) {
				s.merge(slot[0]);
				memcpy(&slot[0], &s, sizeof(Slot));
			}
			else if (s < slot[1]) {
				s.merge(slot[1]);
				memcpy(&slot[1], &s, sizeof(Slot));
			}
			else {
				slot[0].merge(slot[1]);
				memcpy(&slot[1], &s, sizeof(Slot));
			}
			std::sort(slot, slot + cell_num);*/
		}
		
	}
	uint64_t query(double w) {
		int32_t m = f * w;
		// std::cout << f << "\n";
		// for (int i = 0; i < slot_num; ++i) {
		// 	std::cout << "[" << pow(log_base, slot[i].start) << ", " << pow(log_base, slot[i].end) << "]: " << slot[i].f << "\n"; 
		// }
		for (int i = 0; i < slot_num; ++i) {
			m -= slot[i].f;
			if (m < 0) {
				return (slot[i].start + slot[i].end) / 2;
			}
		}
	}
	uint32_t f = 0;
	uint32_t slot_num = 0;
	uint32_t cell_num;
private:
	Slot* slot;
};

template<typename ID_TYPE, typename DATA_TYPE>
class DDSketch {
public:
	DDSketch() {}
	DDSketch(uint32_t memory) {
		uint32_t total_cell = memory * 1024 / sizeof(Slot);
		cell_num = total_cell;
		bucket_num = 1;
		bucket = new Bucket [bucket_num];
		for (int i = 0; i < bucket_num; ++i) {
			bucket[i].init(cell_num);
		}
	}
	~DDSketch() {}
	void insert(ID_TYPE id, DATA_TYPE time) {
		uint32_t index = hash(id, 1) % bucket_num;
		uint32_t T = (log(time) / log(log_base) + 0.5);
		bucket[index].insert(T);
	}
	DATA_TYPE query(ID_TYPE id, double w) {
		uint32_t index = hash(id, 1) % bucket_num;
		return pow(log_base, bucket[index].query(w));
	}
	Bucket* bucket;
	uint32_t bucket_num;
	uint32_t cell_num;
};

#endif