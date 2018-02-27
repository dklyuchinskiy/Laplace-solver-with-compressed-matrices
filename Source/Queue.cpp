#include "definitions.h"
#include "templates.h"
#include "TestSuite.h"

/*********************************
Implementation of Queue by List.
(We are going to use it later 
in PrintRanksInWidth function)
*********************************/

void init(struct my_queue* &q)
{
	q = (struct my_queue*)malloc(sizeof(struct my_queue));
	q->first = NULL;
	q->last = NULL;
}

bool my_empty(struct my_queue* q)
{
	if (q->first == NULL && q->last == NULL) return true;
	else return false;
}

void push(struct my_queue* &q, mnode* node)
{
	qlist *item = (qlist*)malloc(sizeof(qlist)); // просто перекидываем указатель по листу, память под node уже выделена на верхнем уровне
	item->node = node;
	item->next = NULL;

	// указатель first установлен всегда на первый элемент, а last - двигается по списку по мере добавления элементов, но всегда
	// указывает на последний элемент

	// Если очередь пуста
	if (my_empty(q))
	{ 
		q->first = item;
		q->last = q->first;
	}
	else
	{
#if 0
		if (q->first == q->last)
		{
			q->first->next = item;
			q->last = item;
		}
		else
		{
			q->last->next = item;
			q->last = item;
		}
#else
		q->last->next = item; // кладем в указатель последнего елемента новый item
		q->last = item; // т.к. добавился в конец новый элемент, то переносим указатель last на него
#endif
	}
}

// удаляем первый элемент из очереди
void pop(struct my_queue* &q)
{
	qlist* temp; // память уже выделена в insert под node, просто перекидываем указатель
	if (my_empty(q)) return;
	else
	{
		temp = q->first; // сохраняем старый первый номер
		q->first = q->first->next;
		if (q->first == NULL) q->last = NULL; // для удаления последнего элемента

		free(temp); // удаляем только указатель, память под mnode удаляется на другом уровне
	}
}

mnode* front(struct my_queue* q)
{
	mnode* x;
	x = q->first->node;
	return x;
}

void print_queue(struct my_queue* q)
{
	qlist* h;
	if (my_empty(q))
	{
		printf("Queue is empty!\n");
		return;
	}
	else
	{
		printf("Cur queue: ");
		for (h = q->first; h != NULL; h = h->next)
			printf("%d ", h->node->p);
		printf("\n");
	}
}
