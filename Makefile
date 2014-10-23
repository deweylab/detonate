.PHONY: all
all: ref-eval/finished rsem-eval/finished

ref-eval/finished:
	@echo 
	@echo ===========================================
	@echo = Building REF-EVAL and its dependencies. =
	@echo ===========================================
	@echo 
	cd ref-eval && $(MAKE) && touch finished

rsem-eval/finished:
	@echo 
	@echo ============================================
	@echo = Building RSEM-EVAL and its dependencies. =
	@echo ============================================
	@echo 
	cd rsem-eval && $(MAKE) && touch finished
