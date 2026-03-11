"""Context Jimi — the shared state manager.

Jimi is the central nervous system of the Darwin Council.
All agents read from and write to Jimi. He maintains:
  - The genome object (evolving as agents annotate it)
  - The annotation config
  - A message bus for agent communication
  - Agent status tracking
  - Decision logs from council deliberations
  - Full audit trail

Named after Jimi Hendrix because he orchestrates the experience.
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from typing import Any, TYPE_CHECKING

from darwin.models import AnnotationConfig, Genome
from darwin.utils.logging import get_logger

if TYPE_CHECKING:
    from darwin.agents.base import Agent, AgentMessage, AgentStatus

log = get_logger("darwin.context.jimi")


@dataclass
class DecisionRecord:
    """A record of a council decision."""
    topic: str
    votes: dict[str, str]          # agent_name -> vote
    rationale: dict[str, str]      # agent_name -> reasoning
    outcome: str
    timestamp: float = field(default_factory=time.time)


class ContextJimi:
    """Shared context and message bus for the Darwin Agent Council.

    Thread-safe — agents can run concurrently if needed.
    """

    def __init__(self, config: AnnotationConfig) -> None:
        self.config = config
        self._lock = threading.RLock()

        # Core state
        self.genome: Genome | None = None
        self.results: dict[str, Any] = {}       # agent_name -> result dict
        self.errors: dict[str, str] = {}        # agent_name -> error message

        # Agent registry
        self._agents: dict[str, Agent] = {}
        self._agent_status: dict[str, AgentStatus] = {}

        # Message bus
        self._message_queue: list[AgentMessage] = []
        self._message_history: list[AgentMessage] = []

        # Decision log
        self._decisions: list[DecisionRecord] = []

        # Metadata / scratch space agents can use
        self.scratch: dict[str, Any] = {}

        log.info("[bold magenta]Context Jimi[/] initialized")

    # ── Agent Registry ───────────────────────────────────────

    def register_agent(self, agent: Agent) -> None:
        """Register an agent with the context."""
        with self._lock:
            self._agents[agent.name] = agent
            self._agent_status[agent.name] = agent.status
            log.debug(f"  Registered agent: {agent.name} ({agent.role})")

    def get_agent(self, name: str) -> Agent | None:
        with self._lock:
            return self._agents.get(name)

    def update_agent_status(self, name: str, status: AgentStatus) -> None:
        with self._lock:
            self._agent_status[name] = status

    def get_council_status(self) -> dict[str, str]:
        """Get the status of all agents."""
        with self._lock:
            return {name: status.value for name, status in self._agent_status.items()}

    # ── Genome State ─────────────────────────────────────────

    def set_genome(self, genome: Genome) -> None:
        with self._lock:
            self.genome = genome
            log.info(
                f"  Genome loaded: {genome.name} "
                f"({genome.num_contigs} contigs, {genome.total_length:,} bp)"
            )

    def get_genome(self) -> Genome:
        with self._lock:
            if self.genome is None:
                raise RuntimeError("No genome loaded in context")
            return self.genome

    # ── Results Store ────────────────────────────────────────

    def store_result(self, agent_name: str, result: dict[str, Any]) -> None:
        """Store an agent's output in the shared context."""
        with self._lock:
            self.results[agent_name] = result
            log.debug(f"  Stored result from {agent_name}")

    def get_result(self, agent_name: str) -> dict[str, Any] | None:
        with self._lock:
            return self.results.get(agent_name)

    def store_error(self, agent_name: str, error: str) -> None:
        with self._lock:
            self.errors[agent_name] = error

    # ── Message Bus ──────────────────────────────────────────

    def route_message(self, message: AgentMessage) -> None:
        """Route a message to its recipient(s)."""
        with self._lock:
            self._message_history.append(message)

            if message.recipient == "council":
                # Broadcast to all agents except sender
                for name, agent in self._agents.items():
                    if name != message.sender:
                        agent.receive(message)
            else:
                # Direct message
                agent = self._agents.get(message.recipient)
                if agent:
                    agent.receive(message)
                else:
                    log.warning(
                        f"  Message from {message.sender} to unknown agent "
                        f"'{message.recipient}'"
                    )

    # ── Council Decisions ────────────────────────────────────

    def record_decision(
        self,
        topic: str,
        votes: dict[str, str],
        rationale: dict[str, str],
        outcome: str,
    ) -> DecisionRecord:
        """Record a council decision."""
        decision = DecisionRecord(
            topic=topic,
            votes=votes,
            rationale=rationale,
            outcome=outcome,
        )
        with self._lock:
            self._decisions.append(decision)
        log.info(f"  Council decision: [bold]{topic}[/] → {outcome}")
        return decision

    def get_decisions(self) -> list[DecisionRecord]:
        with self._lock:
            return list(self._decisions)

    # ── Audit Trail ──────────────────────────────────────────

    def get_full_audit(self) -> dict:
        """Get the complete audit trail of the annotation run."""
        with self._lock:
            return {
                "config": {
                    "input": str(self.config.input_path),
                    "output": str(self.config.output_dir),
                    "locus_tag": self.config.locus_tag_prefix,
                },
                "agents": {
                    name: {
                        "status": self._agent_status.get(name, "unknown"),
                        "has_result": name in self.results,
                        "has_error": name in self.errors,
                    }
                    for name in self._agents
                },
                "decisions": [
                    {
                        "topic": d.topic,
                        "outcome": d.outcome,
                        "votes": d.votes,
                    }
                    for d in self._decisions
                ],
                "message_count": len(self._message_history),
                "genome_summary": self.genome.summary() if self.genome else None,
            }
